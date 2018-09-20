! ======================================================================
! A basis of states which spans a subspace.
! ======================================================================
module subspace_basis_module
  use common_module
  
  use monomial_state_module
  use braket_module
  implicit none
  
  private
  
  public :: SubspaceBasis
  public :: SubspaceWavevectorBasis
  public :: size
  public :: generate_subspace_basis
  
  ! Conversions between bases of states.
  !type, extends(Stringable) :: StateConversion
  !  integer,  allocatable :: ids(:)
  !  real(dp), allocatable :: coefficients(:)
  !end type
  
  ! The states at a single wavevector.
  type, extends(Stringsable) :: SubspaceWavevectorBasis
    ! The ID and frequency of the subspace.
    integer  :: subspace_id
    real(dp) :: frequency
    ! The wavevector.
    type(FractionVector) :: wavevector
    ! The states at the wavevector.
    type(MonomialState), allocatable :: states(:)
    ! Conversion matrices from the stored states(:) (which are normalised
    !    but not in general orthogonal) to an orthonormal harmonic basis,
    !    and back again.
    type(RealMatrix), private :: states_to_basis_
    type(RealMatrix), private :: basis_to_states_
    ! The occupations of the harmonic basis states.
    ! e.g. the occupation of |0> is zero, and the occupation of |2,4> is 6.
    integer, allocatable :: harmonic_occupations(:)
    
    type(IntArray1D), allocatable, private :: states_to_basis_ids_(:)
    type(RealVector), allocatable, private :: states_to_basis_coefficients_(:)
    type(IntArray1D), allocatable, private :: basis_to_states_ids_(:)
    type(RealVector), allocatable, private :: basis_to_states_coefficients_(:)
  contains
    ! Set the frequency of the basis.
    procedure, public :: set_frequency => set_frequency_SubspaceWavevectorBasis
    
    ! Transform vectors of coefficients.
    procedure, public :: coefficients_states_to_basis
    procedure, public :: coefficients_basis_to_states
    
    ! Transform operator matrices.
    procedure, public :: operator_states_to_basis
    procedure, public :: operator_basis_to_states
    
    ! I/O.
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
    ! Set the frequency of the basis.
    procedure, public :: set_frequency => set_frequency_SubspaceBasis
    
    ! I/O.
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
   & states_to_basis,basis_to_states,harmonic_occupations) result(this)!,                    &
!   & states_to_basis_ids,states_to_basis_coefficients,basis_to_states_ids,    &
!   & basis_to_states_coefficients) result(this)
  implicit none
  
  integer,              intent(in) :: subspace_id
  real(dp),             intent(in) :: frequency
  type(FractionVector), intent(in) :: wavevector
  type(MonomialState),  intent(in) :: states(:)
  type(RealMatrix),     intent(in) :: states_to_basis
  type(RealMatrix),     intent(in) :: basis_to_states
  integer,              intent(in) :: harmonic_occupations(:)
!  type(IntArray1D),     intent(in) :: states_to_basis_ids(:)
!  type(RealVector),     intent(in) :: states_to_basis_coefficients(:)
!  type(IntArray1D),     intent(in) :: basis_to_states_ids(:)
!  type(RealVector),     intent(in) :: basis_to_states_coefficients(:)
  type(SubspaceWavevectorBasis)    :: this
  
  if (any(states%subspace_id/=subspace_id)) then
    call print_line(CODE_ERROR//": Subspaces don't match.")
    call err()
  endif
  
  this%subspace_id                   = subspace_id
  this%frequency                     = frequency
  this%wavevector                    = wavevector
  this%states                        = states
  this%states_to_basis_              = states_to_basis
  this%basis_to_states_              = basis_to_states
  this%harmonic_occupations          = harmonic_occupations
!  this%states_to_basis_ids_          = states_to_basis_ids
!  this%states_to_basis_coefficients_ = states_to_basis_coefficients
!  this%basis_to_states_ids_          = basis_to_states_ids
!  this%basis_to_states_coefficients_ = basis_to_states_coefficients
end function

function new_SubspaceBasis(subspace_id,frequency,wavevectors) result(this)
  implicit none
  
  integer,                       intent(in) :: subspace_id
  real(dp),                      intent(in) :: frequency
  type(SubspaceWavevectorBasis), intent(in) :: wavevectors(:)
  type(SubspaceBasis)                       :: this
  
  this%subspace_id = subspace_id
  this%frequency   = frequency
  this%wavevectors = wavevectors
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

! Set the frequency of the basis.
impure elemental subroutine set_frequency_SubspaceBasis(this,frequency)
  implicit none
  
  class(SubspaceBasis), intent(inout) :: this
  real(dp),             intent(in)    :: frequency
  
  this%frequency = frequency
  call this%wavevectors%set_frequency(frequency)
end subroutine

impure elemental subroutine set_frequency_SubspaceWavevectorBasis(this, &
   & frequency)
  implicit none
  
  class(SubspaceWavevectorBasis), intent(inout) :: this
  real(dp),                       intent(in)    :: frequency
  
  this%frequency = frequency
  this%states%frequency = frequency
end subroutine

! Generates states up to a given power.
function generate_subspace_basis(subspace,frequency,modes,qpoints, &
   & supercell,maximum_power) result(output)
  implicit none
  
  type(DegenerateSubspace), intent(in) :: subspace
  real(dp),                 intent(in) :: frequency
  type(ComplexMode),        intent(in) :: modes(:)
  type(QpointData),         intent(in) :: qpoints(:)
  type(StructureData),      intent(in) :: supercell
  integer,                  intent(in) :: maximum_power
  type(SubspaceBasis)                  :: output
  
  ! States in the basis.
  type(MonomialState), allocatable :: states(:)
  
  ! Variables for matching wavevectors to states.
  type(QpointData), allocatable :: wavevectors(:)
  integer,          allocatable :: sort_key(:)
  type(QpointData), allocatable :: unique_wavevectors(:)
  type(IntArray1D), allocatable :: wavevector_state_ids(:)
  
  ! Variables for orthonormalising states.
  type(MonomialState),       allocatable :: wavevector_states(:)
  real(dp),                  allocatable :: states_to_basis(:,:)
  real(dp),                  allocatable :: basis_to_states(:,:)
  integer,                   allocatable :: harmonic_occupations(:)
  type(MonomialState),       allocatable :: unique_states(:)
  integer,                   allocatable :: matching_state_ids(:)
  type(MonomialState),       allocatable :: matching_states(:)
  real(dp),                  allocatable :: state_overlaps(:,:)
  real(dp),                  allocatable :: harmonic_hamiltonian(:,:)
  type(SymmetricEigenstuff), allocatable :: estuff(:)
  real(dp),                  allocatable :: matching_states_to_basis(:,:)
  real(dp),                  allocatable :: matching_basis_to_states(:,:)
  real(dp),                  allocatable :: matching_basis_to_harmonic(:,:)
  real(dp),                  allocatable :: matching_harmonic_to_basis(:,:)
  integer                                :: state
  
  type(IntArray1D), allocatable :: states_to_basis_ids(:)
  type(IntArray1D), allocatable :: basis_to_states_ids(:)
  type(RealVector), allocatable :: states_to_basis_coefficients(:)
  type(RealVector), allocatable :: basis_to_states_coefficients(:)
  
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
    allocate( states_to_basis( size(wavevector_state_ids(i)),               &
            &           size(wavevector_state_ids(i))  ),                   &
            & basis_to_states( size(wavevector_state_ids(i)),               &
            &            size(wavevector_state_ids(i))  ),                  &
            & harmonic_occupations(size(wavevector_state_ids(i))),          &
            & states_to_basis_ids(size(wavevector_state_ids(i))),           &
            & basis_to_states_ids(size(wavevector_state_ids(i))),           &
            & states_to_basis_coefficients(size(wavevector_state_ids(i))),  &
            & basis_to_states_coefficients(size(wavevector_state_ids(i))),  &
            & stat=ialloc); call err(ialloc)
    states_to_basis = 0
    basis_to_states = 0
    unique_states = wavevector_states(set(wavevector_states, states_overlap))
    state = 0
    ! Loop over the sets of states with non-zero overlap.
    do j=1,size(unique_states)
      ! Select the set of states with non-zero overlap.
      matching_state_ids = filter( wavevector_states, &
                                 & states_overlap,    &
                                 & unique_states(j)   )
      matching_states = wavevector_states(matching_state_ids)
      allocate( state_overlaps( size(matching_states),              &
              &                 size(matching_states)),             &
              & harmonic_hamiltonian( size(matching_states),        &
              &                       size(matching_states)),       &
              & matching_states_to_basis( size(matching_states),    &
              &                           size(matching_states)),   &
              & matching_basis_to_states( size(matching_states),    &
              &                           size(matching_states)),   &
              & matching_basis_to_harmonic( size(matching_states),  &
              &                             size(matching_states)), &
              & matching_harmonic_to_basis( size(matching_states),  &
              &                             size(matching_states)), &
              & stat=ialloc); call err(ialloc)
      
      ! Construct the overlap matrix between the states.
      do k=1,size(matching_states)
        do l=1,size(matching_states)
          state_overlaps(l,k) = braket(matching_states(l), matching_states(k))
          harmonic_hamiltonian(l,k) = harmonic_potential_energy(          &
                                  &          matching_states(l),          &
                                  &          matching_states(k),          &
                                  &          subspace,                    &
                                  &          supercell           )        &
                                  & + kinetic_energy( matching_states(l), &
                                  &                   matching_states(k), &
                                  &                   subspace,           &
                                  &                   supercell           )
        enddo
      enddo
      
      ! Diagonalise the overlap matrix, and construct the mapping between
      !    the monomial states and an orthonormal basis.
      estuff = diagonalise_symmetric(state_overlaps)
      do k=1,size(estuff)
        matching_states_to_basis(k,:) = estuff(k)%evec / sqrt(estuff(k)%eval)
        matching_basis_to_states(:,k) = estuff(k)%evec * sqrt(estuff(k)%eval)
      enddo
      
      ! Transform the harmonic Hamiltonian into this orthonormal basis.
      harmonic_hamiltonian = dble( mat(matching_states_to_basis)            &
                               & * mat(harmonic_hamiltonian)                &
                               & * transpose(mat(matching_states_to_basis)) )
      
      ! Diagonalise the Hamiltonian in the orthonormal basis, and construct
      !    the mapping from the orthonromal basis to the harmonic basis.
      ! N.B. the harmonic basis is also orthonormal.
      estuff = diagonalise_symmetric(harmonic_hamiltonian)
      do k=1,size(estuff)
        matching_basis_to_harmonic(k,:) = estuff(k)%evec
        matching_harmonic_to_basis(:,k) = estuff(k)%evec
      enddo
      
      ! Transform the mapping between the monomial states and the orthonormal
      !    basis so that it now maps between the monomial states and the
      !    harmonic basis.
      matching_states_to_basis = dble( mat(matching_basis_to_harmonic) &
                                   & * mat(matching_states_to_basis)   )
      matching_basis_to_states = dble( mat(matching_basis_to_states)   &
                                   & * mat(matching_harmonic_to_basis) )
      
      ! Fill in output data.
      states_to_basis(matching_state_ids,matching_state_ids) = &
         & matching_states_to_basis
      basis_to_states(matching_state_ids,matching_state_ids) = &
         & matching_basis_to_states
      harmonic_occupations(matching_state_ids) =           &
         & nint( estuff%eval*supercell%sc_size / frequency &
         &     - 0.5_dp*size(subspace)                     )
      
      ! Deallocate temporary arrays.
      deallocate( state_overlaps,             &
                & harmonic_hamiltonian,       &
                & matching_states_to_basis,   &
                & matching_basis_to_states,   &
                & matching_basis_to_harmonic, &
                & matching_harmonic_to_basis, &
                & stat=ialloc); call err(ialloc)
    enddo
    
    do j=1,size(states_to_basis,1)
      states_to_basis_ids(j) = filter(abs(states_to_basis(j,:))>1e-10_dp)
      states_to_basis_coefficients(j) = &
         & states_to_basis(j,states_to_basis_ids(j)%i)
      basis_to_states_ids(j) = filter(abs(basis_to_states(j,:))>1e-10_dp)
      basis_to_states_coefficients(j) = &
         & basis_to_states(j,basis_to_states_ids(j)%i)
    enddo
    
    wavevector_bases(i) = SubspaceWavevectorBasis( &
                   & subspace%id,                  &
                   & frequency,                    &
                   & unique_wavevectors(i)%qpoint, &
                   & wavevector_states,            &
                   & mat(states_to_basis),         &
                   & mat(basis_to_states),         &
                   & harmonic_occupations)!,         &
!                   & states_to_basis_ids,          &
!                   & states_to_basis_coefficients, &
!                   & basis_to_states_ids,          &
!                   & basis_to_states_coefficients  )
    
    deallocate( states_to_basis,              &
              & basis_to_states,              &
              & harmonic_occupations,         &
              & states_to_basis_ids,          &
              & states_to_basis_coefficients, &
              & basis_to_states_ids,          &
              & basis_to_states_coefficients, &
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
    
    select type(bra); type is(MonomialState)
      select type(ket); type is(MonomialState)
        output = finite_overlap(bra,ket)
      end select
    end select
  end function
end function

! ----------------------------------------------------------------------
! Transform coefficients between monomial states and the orthonormal basis.
! ----------------------------------------------------------------------
! N.B. the private matrices states_to_basis_ and basis_to_states_ define how
!    the states themselves transform. The coefficients transform as the inverse
!    transpose of this.
function coefficients_states_to_basis(this,coefficients) result(output)
  implicit none
  
  class(SubspaceWavevectorBasis), intent(in) :: this
  real(dp),                       intent(in) :: coefficients(:)
  real(dp), allocatable                      :: output(:)
  
  output = dble(transpose(this%basis_to_states_) * vec(coefficients))
end function

function coefficients_basis_to_states(this,coefficients) result(output)
  implicit none
  
  class(SubspaceWavevectorBasis), intent(in) :: this
  real(dp),                       intent(in) :: coefficients(:)
  real(dp), allocatable                      :: output(:)
  
  output = dble(transpose(this%states_to_basis_) * vec(coefficients))
end function

! ----------------------------------------------------------------------
! Transform operators between monomial states and the orthonormal basis.
! ----------------------------------------------------------------------
function operator_states_to_basis(this,operator) result(output)
  implicit none
  
  class(SubspaceWavevectorBasis), intent(in) :: this
  real(dp),                       intent(in) :: operator(:,:)
  real(dp), allocatable                      :: output(:,:)
  
  output = dble( this%states_to_basis_            &
             & * mat(operator)                    &
             & * transpose(this%states_to_basis_) )
end function

function operator_basis_to_states(this,operator) result(output)
  implicit none
  
  class(SubspaceWavevectorBasis), intent(in) :: this
  real(dp),                       intent(in) :: operator(:,:)
  real(dp), allocatable                      :: output(:,:)
  
  output = dble( this%basis_to_states_            &
             & * mat(operator)                    &
             & * transpose(this%basis_to_states_) )
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
  type(MonomialState), allocatable :: states(:)
  integer,             allocatable :: harmonic_occupations(:)
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
    
    no_states = (size(input)-7)/4
    
    allocate(states(no_states), stat=ialloc); call err(ialloc)
    do i=1,no_states
      line = split_line(input(4+i))
      states(i) = MonomialState([input(1:2), str('State'), line(3)])
    enddo
    
    allocate(harmonic_occupations(no_states), stat=ialloc); call err(ialloc)
    do i=1,no_states
      line = split_line(input(5+no_states+i))
      harmonic_occupations(i) = int(line(3))
    enddo
    
    states_to_basis = RealMatrix(input(7+2*no_states:6+3*no_states))
    basis_to_states = RealMatrix(input(8+3*no_states:7+4*no_states))
    
    this = SubspaceWavevectorBasis( subspace_id,         &
                                  & frequency,           &
                                  & wavevector,          &
                                  & states,              &
                                  & states_to_basis,     &
                                  & basis_to_states,     &
                                  & harmonic_occupations )
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
      output = [output, '|'//i//'> = '//state_strings(size(state_strings))]
    enddo
    output = [ output,                     &
             & str('Harmonic occupations') ]
    do i=1,size(this%harmonic_occupations)
      output = [ output,                                      &
               & '|'//i//'> : '//this%harmonic_occupations(i) ]
    enddo
    output = [ output,                            &
             & str('States to Basis Conversion'), &
             & str(this%states_to_basis_),        &
             & str('Basis to States Conversion'), &
             & str(this%basis_to_states_)         ]
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
    line = split_line(input(1))
    subspace_id = int(line(2))
    
    line = split_line(input(2))
    frequency = dble(line(2))
    
    starting_lines = [integer::]
    do i=3,size(input)
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
