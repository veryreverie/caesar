! ======================================================================
! A basis of states which spans a subspace.
! ======================================================================
module subspace_basis_module
  use common_module
  
  use monomial_state_module
  use braket_module
  use state_conversion_module
  use wavevector_basis_module
  
  use old_wavevector_basis_module
  implicit none
  
  private
  
  public :: SubspaceBasis
  public :: size
  public :: generate_subspace_basis
  
  ! All states spanning the subspace.
  type, extends(Stringsable) :: SubspaceBasis
    ! The maximum power of the monomial states.
    ! N.B. this is also the maximum occupation of the harmonic basis states.
    integer  :: maximum_power
    ! The ID and frequency of the subspace.
    integer  :: subspace_id
    real(dp) :: frequency
    ! The states, wavevector by wavevector.
    type(WavevectorBasis), allocatable :: wavevectors(:)
  contains
    ! Set the frequency of the basis.
    procedure, public :: set_frequency => set_frequency_SubspaceBasis
    
    ! I/O.
    procedure, public :: read  => read_SubspaceBasis
    procedure, public :: write => write_SubspaceBasis
  end type
  
  interface SubspaceBasis
    module procedure new_SubspaceBasis
    module procedure new_SubspaceBasis_Strings
    module procedure new_SubspaceBasis_StringArray
  end interface
  
  interface size
    module procedure size_SubspaceBasis
  end interface
contains

! Constructors and size functions.
function new_SubspaceBasis(maximum_power,subspace_id,frequency,wavevectors) &
   & result(this)
  implicit none
  
  integer,               intent(in) :: maximum_power
  integer,               intent(in) :: subspace_id
  real(dp),              intent(in) :: frequency
  type(WavevectorBasis), intent(in) :: wavevectors(:)
  type(SubspaceBasis)               :: this
  
  this%maximum_power = maximum_power
  this%subspace_id   = subspace_id
  this%frequency     = frequency
  this%wavevectors   = wavevectors
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
  
  type(WavevectorBasis), allocatable :: wavevectors(:)
  
  type(SubspaceWavevectorBasis), allocatable :: old_wavevectors(:)
  
  integer :: i,j,k
  real(dp), allocatable :: brakets(:,:)
  
  wavevectors = generate_subspace_basis_states( subspace,     &
                                              & frequency,    &
                                              & modes,        &
                                              & qpoints,      &
                                              & maximum_power )
  
  ! TODO: remove below.
!  old_wavevectors = generate_old_basis(subspace,frequency,modes,qpoints,supercell,maximum_power)
!  
!  do i=1,size(wavevectors)
!    call print_line(repeat('=',50))
!    call print_lines(wavevectors(i))
!    call print_line(repeat('.',50))
!    call print_lines(old_wavevectors(first(old_wavevectors%wavevector==wavevectors(i)%wavevector)))
!    call print_line(repeat('-',50))
!    
!    allocate(brakets(size(wavevectors(i)),size(wavevectors(i))))
!    do j=1,size(wavevectors(i))
!      do k=1,size(wavevectors(i))
!        brakets(j,k) = braket(wavevectors(i)%monomial_states(j), wavevectors(i)%monomial_states(k))
!      enddo
!    enddo
!    brakets = wavevectors(i)%operator_states_to_basis(brakets)
!    call print_line('')
!    call print_lines(mat(brakets))
!    deallocate(brakets)
!  enddo
  ! TODO: remove above.
  
  output = SubspaceBasis( maximum_power, &
                        & subspace%id,   &
                        & frequency,     &
                        & wavevectors    )
end function

! ----------------------------------------------------------------------
! I/O.
! ----------------------------------------------------------------------
subroutine read_SubspaceBasis(this,input)
  implicit none
  
  class(SubspaceBasis), intent(out) :: this
  type(String),         intent(in)  :: input(:)
  
  integer                            :: maximum_power
  integer                            :: subspace_id
  real(dp)                           :: frequency
  type(WavevectorBasis), allocatable :: wavevectors(:)
  
  integer, allocatable :: starting_lines(:)
  integer, allocatable :: ending_lines(:)
  
  type(String), allocatable :: line(:)
  type(String), allocatable :: wavevector_lines(:)
  
  integer :: i,ialloc
  
  select type(this); type is(SubspaceBasis)
    line = split_line(input(1))
    maximum_power = int(line(4))
    
    line = split_line(input(2))
    subspace_id = int(line(3))
    
    line = split_line(input(3))
    frequency = dble(line(3))
    
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
      wavevector_lines = [input(:3), input(starting_lines(i):ending_lines(i))]
      wavevectors(i) = WavevectorBasis(wavevector_lines)
    enddo
    
    this = SubspaceBasis(maximum_power,subspace_id,frequency,wavevectors)
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
    output = [ 'Maximum power : '//this%maximum_power, &
             & 'Subspace      : '//this%subspace_id, &
             & 'Frequency     : '//this%frequency   ]
    do i=1,size(this%wavevectors)
      wavevector_strings = str(this%wavevectors(i))
      output = [ output, wavevector_strings(4:) ]
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
