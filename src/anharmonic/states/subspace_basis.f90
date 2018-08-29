! ======================================================================
! A basis of states which spans a subspace.
! ======================================================================
module subspace_basis_module
  use common_module
  
  use subspace_state_module
  implicit none
  
  private
  
  public :: SubspaceBasis
  public :: SubspaceWavevectorBasis
  public :: size
  public :: generate_subspace_basis
  
  ! The states at a single wavevector.
  type, extends(Stringsable) :: SubspaceWavevectorBasis
    ! The ID and frequency of the subspace.
    integer  :: subspace_id
    real(dp) :: frequency
    ! The wavevector.
    type(FractionVector) :: wavevector
    ! The states at the wavevector.
    type(SubspaceState), allocatable :: states(:)
    ! Conversion matrices from the stored states(:) (which are normalised
    !    but not in general orthogonal) to an orthonormal basis,
    !    and back again.
    type(RealMatrix) :: states_to_basis
    type(RealMatrix) :: basis_to_states
  contains
    procedure, public :: read  => read_SubspaceWavevectorBasis
    procedure, public :: write => write_SubspaceWavevectorBasis
  end type
  
  ! All states spanning the subspace.
  type, extends(Stringsable) :: SubspaceBasis
    ! The ID and frequency of the subspace.
    integer  :: subspace_id
    real(dp) :: frequency
    ! The states, wavevector by wavevector.
    type(SubspaceWavevectorBasis), allocatable :: wavevectors(:)
  contains
    procedure, public :: read  => read_SubspaceBasis
    procedure, public :: write => write_SubspaceBasis
  end type
  
  interface SubspaceWavevectorBasis
    module procedure new_SubspaceWavevectorBasis
    module procedure new_SubspaceWavevectorBasis_Strings
    module procedure new_SubspaceWavevectorBasis_StringArray
  end interface
  
  interface SubspaceBasis
    module procedure new_SubspaceBasis
    module procedure new_SubspaceBasis_Strings
    module procedure new_SubspaceBasis_StringArray
  end interface
  
  interface size
    module procedure size_SubspaceWavevectorBasis
    module procedure size_SubspaceBasis
  end interface
contains

! Constructors and size functions.
function new_SubspaceWavevectorBasis(subspace_id,frequency,wavevector,states, &
   & states_to_basis,basis_to_states) result(this)
  implicit none
  
  integer,              intent(in) :: subspace_id
  real(dp),             intent(in) :: frequency
  type(FractionVector), intent(in) :: wavevector
  type(SubspaceState),  intent(in) :: states(:)
  type(RealMatrix),     intent(in) :: states_to_basis
  type(RealMatrix),     intent(in) :: basis_to_states
  type(SubspaceWavevectorBasis)    :: this
  
  if (any(states%subspace_id/=subspace_id)) then
    call print_line(CODE_ERROR//": Subspaces don't match.")
    call err()
  endif
  
  this%subspace_id     = subspace_id
  this%frequency       = frequency
  this%wavevector      = wavevector
  this%states          = states
  this%states_to_basis = states_to_basis
  this%basis_to_states = basis_to_states
end function

function new_SubspaceBasis(subspace_id,frequency,wavevectors) result(this)
  implicit none
  
  integer,                       intent(in) :: subspace_id
  real(dp),                      intent(in) :: frequency
  type(SubspaceWavevectorBasis), intent(in) :: wavevectors(:)
  type(SubspaceBasis)                       :: this
  
  this%subspace_id          = subspace_id
  this%frequency            = frequency
  this%wavevectors          = wavevectors
end function

function size_SubspaceWavevectorBasis(this) result(output)
  implicit none
  
  type(SubspaceWavevectorBasis), intent(in) :: this
  integer                                   :: output
  
  output = size(this%states)
end function

function size_SubspaceBasis(this) result(output)
  implicit none
  
  type(SubspaceBasis), intent(in) :: this
  integer                         :: output
  
  output = size(this%wavevectors)
end function

! Generates states up to a given power.
function generate_subspace_basis(subspace,frequency,modes,qpoints, &
   & maximum_power) result(output)
  implicit none
  
  type(DegenerateSubspace), intent(in) :: subspace
  real(dp),                 intent(in) :: frequency
  type(ComplexMode),        intent(in) :: modes(:)
  type(QpointData),         intent(in) :: qpoints(:)
  integer,                  intent(in) :: maximum_power
  type(SubspaceBasis)                  :: output
  
  ! States in the basis.
  type(SubspaceState), allocatable :: states(:)
  
  ! Variables for matching wavevectors to states.
  type(QpointData), allocatable :: wavevectors(:)
  integer,          allocatable :: sort_key(:)
  type(QpointData), allocatable :: unique_wavevectors(:)
  type(IntArray1D), allocatable :: wavevector_state_ids(:)
  
  ! Variables for orthonormalising states.
  type(SubspaceState),       allocatable :: wavevector_states(:)
  real(dp),                  allocatable :: states_to_basis(:,:)
  real(dp),                  allocatable :: basis_to_states(:,:)
  type(SubspaceState),       allocatable :: unique_states(:)
  integer,                   allocatable :: matching_state_ids(:)
  type(SubspaceState),       allocatable :: matching_states(:)
  real(dp),                  allocatable :: state_overlaps(:,:)
  type(SymmetricEigenstuff), allocatable :: estuff(:)
  integer                                :: state
  
  ! Output variables.
  type(SubspaceWavevectorBasis), allocatable :: wavevector_bases(:)
  
  ! Temporary variables.
  integer :: i,j,k,l,ialloc
  
  ! ------------------------------
  ! Generate all states in the basis.
  ! ------------------------------
  states = generate_subspace_states( subspace,     &
                                   & frequency,    &
                                   & modes,        &
                                   & maximum_power )
  
  ! ------------------------------
  ! Organise states by wavevector.
  ! ------------------------------
  ! Calculate the wavevector of each state.
  wavevectors = [(states(i)%wavevector(modes,qpoints), i=1, size(states))]
  
  ! Sort states and wavevectors by wavevector ID.
  sort_key = sort(wavevectors%id)
  states = states(sort_key)
  wavevectors = wavevectors(sort_key)
  
  ! Generate a de-duplicated list of wavevectors.
  unique_wavevectors = wavevectors(set(wavevectors%id))
  
  ! Generate a list of which states correspond to which wavevector.
  wavevector_state_ids = [(                                     &
     & array(filter(wavevectors%id==unique_wavevectors(i)%id)), &
     & i=1,                                                     &
     & size(unique_wavevectors)                                 )]
  
  ! ------------------------------
  ! Generate conversions to and from orthonormal basis from states.
  ! ------------------------------
  allocate( wavevector_bases(size(unique_wavevectors)), &
          & stat=ialloc); call err(ialloc)
  do i=1,size(unique_wavevectors)
    wavevector_states = states(wavevector_state_ids(i)%i)
    allocate( states_to_basis( size(wavevector_state_ids(i)), &
            &           size(wavevector_state_ids(i))  ),     &
            & basis_to_states( size(wavevector_state_ids(i)), &
            &            size(wavevector_state_ids(i))  ),    &
            & stat=ialloc); call err(ialloc)
    states_to_basis = 0
    basis_to_states = 0
    unique_states = wavevector_states(set(wavevector_states, states_overlap))
    state = 0
    do j=1,size(unique_states)
      matching_state_ids = filter( wavevector_states, &
                                 & states_overlap,    &
                                 & unique_states(j)   )
      matching_states = wavevector_states(matching_state_ids)
      allocate( state_overlaps(size(matching_states), size(matching_states)), &
              & stat=ialloc); call err(ialloc)
      do k=1,size(matching_states)
        do l=1,size(matching_states)
          state_overlaps(l,k) = braket(matching_states(l), matching_states(k))
        enddo
      enddo
      estuff = diagonalise_symmetric(state_overlaps)
      do k=1,size(estuff)
        state = state+1
        states_to_basis(state, matching_state_ids) = estuff(k)%evec &
                                                 & / sqrt(estuff(k)%eval)
        basis_to_states(matching_state_ids, state) = estuff(k)%evec &
                                                 & * sqrt(estuff(k)%eval)
      enddo
      deallocate(state_overlaps)
    enddo
    
    wavevector_bases(i) = SubspaceWavevectorBasis( &
                   & subspace%id,                  &
                   & frequency,                    &
                   & unique_wavevectors(i)%qpoint, &
                   & wavevector_states,            &
                   & mat(states_to_basis),         &
                   & mat(basis_to_states)          )
    
    deallocate( states_to_basis, &
              & basis_to_states, &
              & stat=ialloc); call err(ialloc)
  enddo
  
  ! ------------------------------
  ! Generate output.
  ! ------------------------------
  output = SubspaceBasis( subspace%id,     &
                        & frequency,       &
                        & wavevector_bases )
contains
  ! Lambda for testing if two states have non-zero overlap.
  function states_overlap(bra,ket) result(output)
    implicit none
    
    class(*), intent(in) :: bra
    class(*), intent(in) :: ket
    logical              :: output
    
    select type(bra); type is(SubspaceState)
      select type(ket); type is(SubspaceState)
        output = finite_overlap(bra,ket)
      end select
    end select
  end function
end function

! ----------------------------------------------------------------------
! I/O.
! ----------------------------------------------------------------------
subroutine read_SubspaceWavevectorBasis(this,input)
  implicit none
  
  class(SubspaceWavevectorBasis), intent(out) :: this
  type(String),                   intent(in)  :: input(:)
  
  integer                          :: subspace_id
  real(dp)                         :: frequency
  type(FractionVector)             :: wavevector
  type(SubspaceState), allocatable :: states(:)
  type(RealMatrix)                 :: states_to_basis
  type(RealMatrix)                 :: basis_to_states
  
  type(String), allocatable :: line(:)
  
  integer :: no_states
  
  integer :: i,ialloc
  
  select type(this); type is(SubspaceWavevectorBasis)
    line = split_line(input(1))
    subspace_id = int(line(2))
    
    line = split_line(input(2))
    frequency = dble(line(2))
    
    line = split_line(input(3))
    wavevector = FractionVector(join(line(2:4)))
    
    no_states = (size(input)-6)/3
    
    allocate(states(no_states), stat=ialloc); call err(ialloc)
    do i=1,no_states
      states(i) = SubspaceState([input(1:2), str('State'), input(4+i)])
    enddo
    
    states_to_basis = RealMatrix(input(6+no_states:5+2*no_states))
    basis_to_states = RealMatrix(input(7+2*no_states:6+3*no_states))
    
    this = SubspaceWavevectorBasis( subspace_id,     &
                                  & frequency,       &
                                  & wavevector,      &
                                  & states,          &
                                  & states_to_basis, &
                                  & basis_to_states  )
  class default
    call err()
  end select
end subroutine

function write_SubspaceWavevectorBasis(this) result(output)
  implicit none
  
  class(SubspaceWavevectorBasis), intent(in) :: this
  type(String), allocatable                  :: output(:)
  
  type(String), allocatable :: state_strings(:)
  
  integer :: i
  
  select type(this); type is(SubspaceWavevectorBasis)
    output = [ 'Subspace '//this%subspace_id,  &
             & 'Frequency '//this%frequency,   &
             & 'Wavevector '//this%wavevector, &
             & str('States')                   ]
    do i=1,size(this%states)
      state_strings = str(this%states(i))
      output = [output, state_strings(size(state_strings))]
    enddo
    output = [ output,                            &
             & str('States to Basis Conversion'), &
             & str(this%states_to_basis),         &
             & str('Basis to States Conversion'), &
             & str(this%basis_to_states)          ]
  class default
    call err()
  end select
end function

function new_SubspaceWavevectorBasis_Strings(input) result(this)
  implicit none
  
  type(String), intent(in)      :: input(:)
  type(SubspaceWavevectorBasis) :: this
  
  call this%read(input)
end function

impure elemental function new_SubspaceWavevectorBasis_StringArray(input) &
   & result(this)
  implicit none
  
  type(StringArray), intent(in) :: input
  type(SubspaceWavevectorBasis) :: this
  
  this = SubspaceWavevectorBasis(str(input))
end function

subroutine read_SubspaceBasis(this,input)
  implicit none
  
  class(SubspaceBasis), intent(out) :: this
  type(String),         intent(in)  :: input(:)
  
  integer                                    :: subspace_id
  real(dp)                                   :: frequency
  type(SubspaceWavevectorBasis), allocatable :: wavevectors(:)
  
  integer, allocatable :: starting_lines(:)
  integer, allocatable :: ending_lines(:)
  
  type(String), allocatable :: line(:)
  type(String), allocatable :: wavevector_lines(:)
  
  integer :: i,ialloc
  
  select type(this); type is(SubspaceBasis)
    starting_lines = [integer::]
    do i=1,size(input)
      line = split_line(input(i))
      if (size(line)>0) then
        if (line(1)=='Wavevector') then
          starting_lines = [starting_lines, i]
        endif
      endif
    enddo
    
    ending_lines = [starting_lines(2:)-1, size(input)]
    
    allocate(wavevectors(size(starting_lines)), stat=ialloc); call err(ialloc)
    do i=1,size(starting_lines)
      wavevector_lines = [input(:2), input(starting_lines(i):ending_lines(i))]
      wavevectors(i) = SubspaceWavevectorBasis(wavevector_lines)
    enddo
    
    this = SubspaceBasis(subspace_id,frequency,wavevectors)
  class default
    call err()
  end select
end subroutine

function write_SubspaceBasis(this) result(output)
  implicit none
  
  class(SubspaceBasis), intent(in) :: this
  type(String), allocatable        :: output(:)
  
  type(String), allocatable :: wavevector_strings(:)
  
  integer :: i
  
  select type(this); type is(SubspaceBasis)
    output = [ 'Subspace '//this%subspace_id, &
             & 'Frequency '//this%frequency   ]
    do i=1,size(this%wavevectors)
      wavevector_strings = str(this%wavevectors(i))
      output = [ output, wavevector_strings(3:) ]
    enddo
  class default
    call err()
  end select
end function

function new_SubspaceBasis_Strings(input) result(this)
  implicit none
  
  type(String), intent(in) :: input(:)
  type(SubspaceBasis)      :: this
  
  call this%read(input)
end function

impure elemental function new_SubspaceBasis_StringArray(input) result(this)
  implicit none
  
  type(StringArray), intent(in) :: input
  type(SubspaceBasis)           :: this
  
  this = SubspaceBasis(str(input))
end function
end module
