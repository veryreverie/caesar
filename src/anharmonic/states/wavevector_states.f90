! ======================================================================
! A set of states in a wavevector basis.
! See wavevector_basis.f90 for more details.
! ======================================================================
module wavevector_states_module
  use common_module
  
  use anharmonic_common_module
  
  use wavevector_state_module
  implicit none
  
  private
  
  public :: startup_wavevector_states
  
  public :: WavevectorStates
  
  type, extends(BasisStates) :: WavevectorStates
    type(WavevectorState), allocatable          :: states(:)
    real(dp),              allocatable          :: energies(:)
    real(dp),              allocatable, private :: weights_(:)
    type(FractionVector),  allocatable, private :: wavevectors_(:)
    type(RealMatrix),      allocatable, private :: density_matrices_(:)
  contains
    procedure, public, nopass :: representation => &
                               & representation_WavevectorStates
    
    procedure, public :: density_matrix
    
    ! I/O.
    procedure, public :: read  => read_WavevectorStates
    procedure, public :: write => write_WavevectorStates
  end type
  
  interface WavevectorStates
    module procedure new_WavevectorStates
    module procedure new_WavevectorStates_BasisStates
    module procedure new_WavevectorStates_Strings
    module procedure new_WavevectorStates_StringArray
  end interface
contains

! Startup procedure.
subroutine startup_wavevector_states()
  implicit none
  
  type(WavevectorStates) :: states
  
  call states%startup()
end subroutine

! ----------------------------------------------------------------------
! WavevectorStates methods.
! ----------------------------------------------------------------------
! Constructors.
function new_WavevectorStates(subspace_id,states,energies,weights) &
   & result(this)
  implicit none
  
  integer,               intent(in)           :: subspace_id
  type(WavevectorState), intent(in)           :: states(:)
  real(dp),              intent(in)           :: energies(:)
  real(dp),              intent(in), optional :: weights(:)
  type(WavevectorStates)                      :: this
  
  integer, allocatable :: wavevector_set(:)
  integer, allocatable :: wavevector_states(:)
  integer              :: no_states
  real(dp), allocatable :: matrix(:,:)
  
  integer :: i,j,k,ialloc
  
  this%subspace_id = subspace_id
  this%states = states
  this%energies = energies
  
  if (present(weights)) then
    ! TODO: Try including sparseness.
    this%weights_ = weights
    
    ! Calculate density matrices.
    wavevector_set = set(states%wavevector, compare_wavevectors)
    this%wavevectors_ = states(wavevector_set)%wavevector
    
    allocate( this%density_matrices_(size(this%wavevectors_)), &
            & stat=ialloc); call err(ialloc)
    do i=1,size(this%wavevectors_)
      no_states = size(this%states(wavevector_set(i))%coefficients)
      allocate(matrix(no_states,no_states), stat=ialloc); call err(ialloc)
      matrix = 0
      do j=1,size(states)
        ! Calculate the elements with l<=k.
        if (states(j)%wavevector==this%wavevectors_(i)) then
          do k=1,no_states
            matrix(:,k) = matrix(:,k) &
                      & + states(j)%coefficients(:) &
                      & * states(j)%coefficients(k) &
                      & * weights(j)
          enddo
        endif
      enddo
      
      this%density_matrices_(i) = mat(matrix)
      deallocate(matrix, stat=ialloc); call err(ialloc)
    enddo
  endif
contains
  ! Helper function to compare wavevectors.
  function compare_wavevectors(this,that) result(output)
    implicit none
    
    class(*), intent(in) :: this
    class(*), intent(in) :: that
    logical              :: output
    
    select type(this); type is(FractionVector)
      select type(that); type is(FractionVector)
        output = this==that
      end select
    end select
  end function
end function

recursive function new_WavevectorStates_BasisStates(input) result(this)
  implicit none
  
  class(BasisStates), intent(in) :: input
  type(WavevectorStates)         :: this
  
  select type(input); type is(WavevectorStates)
    this = input
  type is(BasisStatesPointer)
    ! WORKAROUND: ifort doesn't recognise the interface to this function
    !    from within this function, so the full name is used instead.
    !this = WavevectorStates(input%states())
    this = new_WavevectorStates_BasisStates(input%states())
  class default
    call err()
  end select
end function

! Type representation.
impure elemental function representation_WavevectorStates() result(output)
  implicit none
  
  type(String) :: output
  
  output = 'wavevector state'
end function

! The density matrix, in the harmonic basis, at a given wavevector.
impure elemental function density_matrix(this,wavevector) result(output)
  implicit none
  
  class(WavevectorStates), intent(in) :: this
  type(FractionVector),    intent(in) :: wavevector
  type(RealMatrix)                    :: output
  
  integer :: i
  
  if (.not. allocated(this%wavevectors_)) then
    call print_line(CODE_ERROR//': density matrices not calculated.')
    call err()
  endif
  
  i = first(this%wavevectors_==wavevector, default=0)
  if (i==0) then
    output = dblemat(zeroes(0,0))
  else
    output = this%density_matrices_(i)
  endif
end function

! ----------------------------------------------------------------------
! I/O.
! ----------------------------------------------------------------------
subroutine read_WavevectorStates(this,input)
  implicit none
  
  class(WavevectorStates), intent(out) :: this
  type(String),            intent(in)  :: input(:)
  
  type(StringArray), allocatable :: sections(:)
  
  logical :: weights_present
  
  integer                            :: subspace_id
  type(WavevectorState), allocatable :: states(:)
  real(dp),              allocatable :: energies(:)
  real(dp),              allocatable :: weights(:)
  
  integer :: i,ialloc
  
  select type(this); type is(WavevectorStates)
    subspace_id = int(token(input(1),2))
    
    sections = split_into_sections(input(2:), separating_line=repeat('=',50))
    
    if (size(sections)>0) then
      weights_present = token(sections(1)%strings(size(sections(1))), 1) &
                   & == 'Weight'
    endif
    
    allocate( states(size(sections)),   &
            & energies(size(sections)), &
            & stat=ialloc); call err(ialloc)
    if (weights_present) then
      allocate(weights(size(sections)), stat=ialloc); call err(ialloc)
    endif
    do i=1,size(sections)
      if (weights_present) then
        states(i) = WavevectorState(sections(i)%strings(:size(sections(i))-2))
        
        energies(i) = dble(token(sections(i)%strings(size(sections(i))-1), 3))
        weights(i) = dble(token(sections(i)%strings(size(sections(i))), 3))
      else
        states(i) = WavevectorState(sections(i)%strings(:size(sections(i))-1))
        
        energies(i) = dble(token(sections(i)%strings(size(sections(i))), 3))
      endif
    enddo
    this = WavevectorStates(subspace_id, states, energies, weights)
  class default
    call err()
  end select
end subroutine

function write_WavevectorStates(this) result(output)
  implicit none
  
  class(WavevectorStates), intent(in) :: this
  type(String), allocatable           :: output(:)
  
  type(StringArray), allocatable :: sections(:)
  
  integer :: i
  
  select type(this); type is(WavevectorStates)
    if (allocated(this%weights_)) then
      sections = [( StringArray([ str(this%states(i)),               &
                  &               'Energy : '//this%energies(i),     &
                  &               'Weight : '//this%weights_(i)  ]), &
                  & i=1,                                             &
                  & size(this%states)                                )]
    else
      sections = [( StringArray([ str(this%states(i)),              &
                  &               'Energy : '//this%energies(i) ]), &
                  & i=1,                                            &
                  & size(this%states)                               )]
    endif
    output = [ str('Subspace: '//this%subspace_id),          &
             & str(sections, separating_line=repeat('=',50)) ]
  class default
    call err()
  end select
end function

function new_WavevectorStates_Strings(input) result(this)
  implicit none
  
  type(String), intent(in) :: input(:)
  type(WavevectorStates)   :: this
  
  call this%read(input)
end function

impure elemental function new_WavevectorStates_StringArray(input) result(this)
  implicit none
  
  type(StringArray), intent(in) :: input
  type(WavevectorStates)        :: this
  
  this = WavevectorStates(str(input))
end function
end module
